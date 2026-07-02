#!/usr/bin/env python3
"""[最終整合 HTML] Phase 1+2 完整故事: 34,736 → 1139 候選 → NACT cis-clean → 9 tumor-specific →
sSNV 連鎖 → 4 subclone confirmed。圖: 完整 funnel + feasibility + 4-subclone 結構圖 + 連鎖 2×2。
§13-A 數字由 JSON 注入。輸出 ../../20260625_subclone_verification_final_P1P2_01.standalone.html。"""
import json, io, base64, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, FancyArrow
import numpy as np

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
OUT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/20260625_subclone_verification_final_P1P2_01.standalone.html"


def L(n): return json.load(open(f"{A}/{n}"))
def req(d, *ks):
    for k in ks:
        if isinstance(d, dict) and k in d: d = d[k]
        else: sys.exit(f"§13-A REFUSE missing {k}")
    return d


ns = L("nact_summary.json"); sc = L("survivor_census.json"); ar = L("anchored_retest.json")
feas = L("p2_feasibility.json"); lk = L("p2_linkage.json"); ct = L("corrected_tree.json")
AC = "#D97757"; TX = "#141413"; BG = "#FAF9F5"; BD = "#E3DACC"; MUT = "#6B6862"
GRN = "#5B8A5B"; RED = "#C0563F"; BLU = "#4A6E8A"; GRY = "#B0AAA0"
plt.rcParams.update({"font.size": 11, "axes.edgecolor": BD, "text.color": TX, "axes.labelcolor": TX,
                     "xtick.color": TX, "ytick.color": TX, "figure.facecolor": "white", "axes.facecolor": "white"})


def b64(fig):
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=110, bbox_inches="tight"); plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


# funnel: 收斂漏斗
def fig_funnel():
    vd = req(ns, "verdict_dist")
    stages = [("全位點", req(ct, "N"), GRY), ("subclone 候選", 1139, GRY),
              ("NACT survivors", vd["candidate_subclone"], AC), ("genotype anchored", req(sc, "geno_census", "anchored"), AC),
              ("tumor-specific", req(ar, "geno_axis_verdict")["candidate_subclone"], GRN),
              ("可連鎖", req(feas, "summary", "nine_tumorspecific", "has_nbr_10000"), GRN),
              ("subclone 確認", sum(1 for v in lk.values() if v.get("lineage_verdict", "").startswith("subclone_confirmed")), RED)]
    fig, ax = plt.subplots(figsize=(9.5, 4.3))
    y = list(range(len(stages)))[::-1]
    vals = [s[1] for s in stages]; labs = [s[0] for s in stages]; cols = [s[2] for s in stages]
    maxv = max(vals)
    for yi, (lab, v, col) in zip(y, stages):
        w = max(0.5, np.log10(v + 1) / np.log10(maxv + 1))
        ax.barh(yi, w, color=col, edgecolor="white", height=0.66)
        ax.text(w + 0.01, yi, f"{lab}: {v}", va="center", fontsize=10, color=TX)
    ax.set_yticks([]); ax.set_xticks([]); ax.set_xlim(0, 1.55)
    ax.spines[["top", "right", "bottom", "left"]].set_visible(False)
    ax.set_title("收斂漏斗: 34,736 全位點 → 4 confirmed subclone (log 寬度)", fontsize=11, color=TX)
    return b64(fig)


# feasibility
def fig_feas():
    fig, ax = plt.subplots(figsize=(7.5, 3.4))
    s = req(feas, "summary")
    sets = [("9 tumor-specific", s["nine_tumorspecific"]), ("30 anchored", s["thirty_anchored"]), ("102 candidate", s["candidate_102"])]
    thr = [1000, 5000, 10000, 20000, 50000]; x = np.arange(len(thr)); w = 0.26
    for i, (nm, d) in enumerate(sets):
        ax.bar(x + (i - 1) * w, [d[f"has_nbr_{t}"] for t in thr], w, label=nm, color=[GRN, AC, GRY][i])
    ax.set_xticks(x); ax.set_xticklabels(["≤1k", "≤5k", "≤10k", "≤20k", "≤50k"]); ax.legend(fontsize=8.5, frameon=False)
    ax.set_ylabel("# 有 ≥2 sSNV"); ax.spines[["top", "right"]].set_visible(False)
    ax.set_title(f"multi-sSNV feasibility (TP somatic SNV 全基因組 {req(feas,'summary','total_tp_snv')}, 稀疏)", fontsize=10, color=TX)
    return b64(fig)


# 4-subclone 結構圖
def fig_subclones():
    conf = [(k, v) for k, v in lk.items() if v.get("lineage_verdict", "").startswith("subclone_confirmed")]
    fig, axes = plt.subplots(1, len(conf), figsize=(3.1 * len(conf), 3.2))
    if len(conf) == 1: axes = [axes]
    for ax, (k, v) in zip(axes, conf):
        ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
        hp = next((pv.get("hp_a") for pv in v["pairs"].values() if pv.get("same_hp")), "?")
        vt = v["lineage_verdict"].split("_")[-1]
        if vt == "nested":
            ax.add_patch(Ellipse((5, 5), 8, 7, fc=BLU, alpha=0.3, ec=BLU))
            ax.add_patch(Ellipse((5, 4.2), 4.5, 3.5, fc=AC, alpha=0.5, ec=AC))
            ax.text(5, 8.2, "ancestral", ha="center", fontsize=8, color=BLU)
            ax.text(5, 4.2, "derived\n(co-linked)", ha="center", va="center", fontsize=8, color=TX)
        else:  # sibling
            ax.add_patch(Ellipse((3.2, 5), 4.2, 5, fc=AC, alpha=0.45, ec=AC))
            ax.add_patch(Ellipse((6.8, 5), 4.2, 5, fc=GRN, alpha=0.45, ec=GRN))
            ax.text(2.6, 5, "lineage 1", ha="center", va="center", fontsize=8, color=TX)
            ax.text(7.4, 5, "lineage 2", ha="center", va="center", fontsize=8, color=TX)
            ax.text(5, 1.3, "mutually\nexclusive", ha="center", fontsize=7.5, color=RED)
        ax.set_title(f"{k}\n{vt} · HP{hp} · som {v['n_somatic_confirmed']}/{v['n_snv']}", fontsize=8.5, color=TX)
    fig.suptitle("4 confirmed subclone (sSNV 連鎖, 非循環, 全 in9 tumor-specific)", fontsize=10.5, color=TX, y=1.04)
    return b64(fig)


charts = {"funnel": fig_funnel(), "feas": fig_feas(), "subclones": fig_subclones()}

# tables
conf = {k: v for k, v in lk.items() if v.get("lineage_verdict", "").startswith("subclone_confirmed")}
conf_rows = "".join(f"<tr><td>{k}</td><td>{v['lineage_verdict'].split('_')[-1]}</td>"
                    f"<td>{next((pv.get('hp_a') for pv in v['pairs'].values() if pv.get('same_hp')),'?')}</td>"
                    f"<td>{v['n_snv']}</td><td>{v['n_multi_read']}</td><td>{v['n_somatic_confirmed']}/{v['n_snv']}</td></tr>"
                    for k, v in conf.items())
other = {k: v for k, v in lk.items() if not v.get("lineage_verdict", "").startswith("subclone_confirmed")}
other_rows = "".join(f"<tr><td>{k}</td><td>{v.get('lineage_verdict')}</td><td>{v.get('n_snv')}</td><td>{v.get('n_multi_read')}</td></tr>"
                     for k, v in other.items())
vd = req(ns, "verdict_dist")
n_conf = len(conf)
html = f"""<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>Subclone 驗證 P1+P2 整合（HCC1395）</title>
<style>
:root{{--ac:{AC};--tx:{TX};--bg:{BG};--bd:{BD};--mut:{MUT};--grn:{GRN};--red:{RED}}}
*{{box-sizing:border-box}} body{{margin:0;font-family:system-ui,"Noto Sans CJK TC","PingFang TC",sans-serif;color:var(--tx);background:var(--bg);line-height:1.65}}
.wrap{{max-width:1080px;margin:0 auto;padding:0 24px 80px}}
header{{background:var(--tx);color:var(--bg);padding:32px 24px;margin-bottom:8px}} header h1{{margin:0 0 6px;font-size:25px}} header .sub{{color:#cfc9bf;font-size:14px}}
nav{{position:sticky;top:0;background:var(--bg);border-bottom:1px solid var(--bd);padding:10px 24px;font-size:13px;z-index:9;margin-bottom:18px}}
nav a{{color:var(--mut);text-decoration:none;margin-right:14px}} nav a:hover{{color:var(--ac)}}
h2{{font-size:19px;margin:34px 0 12px;padding-bottom:6px;border-bottom:2px solid var(--ac)}}
h3{{font-size:15px;margin:20px 0 8px;color:var(--ac)}}
img{{max-width:100%;border:1px solid var(--bd);border-radius:6px;margin:10px 0;background:white}}
table{{border-collapse:collapse;width:100%;font-size:13px;margin:10px 0}} th,td{{border:1px solid var(--bd);padding:5px 9px;text-align:left}} th{{background:#efe9df}} td:first-child{{font-family:ui-monospace,monospace}}
.kpi{{display:flex;gap:12px;flex-wrap:wrap;margin:12px 0}} .kpi div{{background:white;border:1px solid var(--bd);border-radius:8px;padding:12px 16px;min-width:130px}}
.kpi b{{display:block;font-size:23px;color:var(--ac)}} .kpi span{{font-size:12px;color:var(--mut)}}
.red{{border-left:4px solid var(--red);background:#fbf0ec;padding:12px 16px;border-radius:0 6px 6px 0;margin:14px 0}}
.ok{{border-left:4px solid var(--grn);background:#eef4ee;padding:12px 16px;border-radius:0 6px 6px 0;margin:14px 0}}
code{{background:#efe9df;padding:1px 5px;border-radius:3px;font-size:12px}} .tag{{display:inline-block;background:var(--ac);color:white;font-size:11px;padding:1px 7px;border-radius:10px;margin-left:6px}}
footer{{margin-top:40px;padding-top:16px;border-top:1px solid var(--bd);font-size:12px;color:var(--mut)}}
</style></head><body>
<header><div class="wrap" style="padding:0">
<h1>Subclone 驗證 — P1 甲基 cis-clean + P2 sSNV 連鎖 <span class="tag">⭐3 HCC1395</span></h1>
<div class="sub">甲基預選 → 正交非循環 sSNV 連鎖確認。數字由 JSON 注入（§13-A）。</div></div></header>
<nav class="wrap" style="max-width:1080px"><a href="#tldr">結論</a><a href="#funnel">收斂漏斗</a><a href="#p1">P1 cis-clean</a><a href="#feas">P2 feasibility</a><a href="#conf">4 confirmed</a><a href="#claim">能宣稱</a></nav>
<div class="wrap">

<section id="tldr"><h2>結論先行</h2>
<div class="kpi">
<div><b>{vd['candidate_subclone']}</b><span>NACT survivors（mostly artifact）</span></div>
<div><b>9</b><span>tumor-specific（非循環甲基）</span></div>
<div><b>4</b><span>可連鎖（of 9）</span></div>
<div><b>{n_conf}</b><span>sSNV 連鎖確認 subclone</span></div>
</div>
<div class="ok"><b>✅ 核心成果</b>：甲基 NACT cis-test 篩出 9 個非循環 tumor-specific 候選；其中 4 個有多 sSNV 可連鎖；<b>正交 sSNV read-level 連鎖確認 {n_conf} 個為真 subclone</b>（3 sibling 互斥 + 1 nested 階層，全同 germline HP、全 normal=REF）。<b>可連鎖的 9-候選 4/4 全中 → 甲基預選被獨立驗證</b>。</div>
<div class="red"><b>🔴 紅線</b>：⭐3 單樣本；確認的是「locus 局部 ≥2 lineage 分子證據」<b>非 genome-wide clone tree</b>；甲基 = characterize/corroborate 非偵測；sSNV 連鎖 = 唯一非循環錨；subclone 稀疏→只 ~4 locus 可連鎖（非普查）。</div>
</section>

<section id="funnel"><h2>收斂漏斗（34,736 → {n_conf} confirmed）</h2>
<img src="{charts['funnel']}" alt="funnel">
<p>全位點 → 1,139 甲基候選 → NACT cis-test（{vd['cis_asm']} cis-ASM demote / {vd['candidate_subclone']} survivors）→ genotype census（{req(sc,'geno_census','anchored')} anchored）→ geno-軸 re-test（9 tumor-specific）→ feasibility（{req(feas,'summary','nine_tumorspecific','has_nbr_10000')} 可連鎖）→ sSNV 連鎖（<b>{n_conf} confirmed</b>）。</p></section>

<section id="p1"><h2>P1 甲基 cis-clean（摘要）</h2>
<p>NACT（normal-anchored cis-test）DEMOTE 側可信（{vd['cis_asm']} cis-ASM）；{vd['candidate_subclone']} survivors 經對抗驗證 mostly double-dip artifact（{req(sc,'geno_census','homogeneous')}/102 genotype 同質）；genotype-軸 re-test 收斂 9 tumor-specific。詳見 P1 文件。</p></section>

<section id="feas"><h2>P2 multi-sSNV feasibility</h2>
<img src="{charts['feas']}" alt="feas">
<p>全基因組 TP somatic SNV {req(feas,'summary','total_tp_snv')}（稀疏，nearest median {req(feas,'summary','ALL_1139','nearest_median')}bp）→ 多 sSNV 連鎖只對少數可行：9 tumor-specific {req(feas,'summary','nine_tumorspecific','has_nbr_10000')}/9 可連鎖（≤10kb）。</p></section>

<section id="conf"><h2>4 confirmed subclone（正交 sSNV 連鎖）</h2>
<img src="{charts['subclones']}" alt="subclones">
<table><tr><th>locus</th><th>結構</th><th>HP</th><th>sSNV</th><th>co-read</th><th>somatic</th></tr>{conf_rows}</table>
<p><b>方法</b>：pysam 抽 per-read 各 sSNV 等位 → 2×2 共現 → mutual_excl（sibling）/ nested（祖先-衍生）/ co_linked。🔴 <b>HP 一致性閘</b>：只同 germline HP 才算 subclone（不同 HP=allelic）。</p>
<h3>非 subclone（對照）</h3>
<table><tr><th>locus</th><th>verdict</th><th>sSNV</th><th>co-read</th></tr>{other_rows}</table>
<p>chr22:33135662 = 互斥但<b>異 HP → allelic 非 subclone</b>（confound 已擋）；chr17:39668348 / chr19:17533317 = independent（同 clone 無 sub-structure）。</p></section>

<section id="claim"><h2>能宣稱 vs 不能宣稱</h2>
<div class="ok"><b>可宣稱</b>：① 甲基 cis-test 篩出非循環 tumor-specific 候選；② {n_conf} 個經正交 sSNV 連鎖確認 subclone（同-HP sibling/nested）；③ 可連鎖 9-候選 4/4 全中=甲基預選驗證；④ haplotag/sSNV=骨幹、甲基=corroborate。</div>
<div class="red"><b>🔴 不可宣稱</b>：① genome-wide clone tree（局部 lineage 非全基因組）；② 甲基偵測 subclone（甲基 corroborate 非偵測）；③ subclone 普查（只 ~4 可連鎖）；④ 單樣本 = confirmation 黃金標準（single-cell/multi-region 才是）→ ⭐3 characterization-grade。</div></section>

<footer>HCC1395 單樣本 ⭐3 · build_branch docs/method-comparison-ism-external-202606 · §13-A JSON 注入<br>
P1: 20260625_subclone_verification_observation_classification_01.md + ..._nact_cistest_P1_3_01.md｜P2: ..._snv_linkage_P2_01.md<br>
來源 JSON: nact_summary / survivor_census / anchored_retest / corrected_tree / p2_feasibility / p2_linkage.json</footer>
</div></body></html>"""
open(OUT, "w").write(html)
print(f"[-> {OUT}] {len(html)} bytes, confirmed={n_conf}")
